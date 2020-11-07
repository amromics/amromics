/* eslint-disable */
import Vue from 'vue'
import Router from 'vue-router'
import SingleSample from '@/views/SingleSample'
import Collection from '@/views/Collection'
import CollectionList from '@/views/CollectionList'
Vue.use(Router)

export default new Router({
  mode: 'history',
  routes: [{ 
      path: '/:cid/:sid',
      component: SingleSample

    },
    {
      path: '/:cid',
      component: Collection
    },
    {
      path: '/',
      component: CollectionList
    }
  ],
})
